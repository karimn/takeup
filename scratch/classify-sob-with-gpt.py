import os
import openai
import pandas as pd
from tenacity import (
    retry,
    stop_after_attempt,
    wait_random_exponential
)
openai.api_key_path=".openai-key"

@retry(wait=wait_random_exponential(min = 1, max = 60), stop = stop_after_attempt(5))
def completion_with_backoff(**kwargs):
    return openai.ChatCompletion.create(**kwargs)

def classify_reason_with_gpt(reason, write_csv = False):
    completion = completion_with_backoff(
        model="gpt-3.5-turbo",
        messages=[
            {
                "role": "user", 
                "content": """
                You will be passed a string explaining why someone may or may not know I was dewormed.
                Classify the reason into the following categories, examples are provided after each category. Only 
                return the category name, not the example. If none of the categories fit well, return 'other'.:

                'campaign' - There was a large deworming campaign in the area, which meant many people went at the same time. e.g.: 'didnt see me there', 'it was announced', 'found me at the treatment point'
                'communication' - I informed the other person. e.g.: 'I told them I was dewormed', 'he informed him'
                'relationship' - I do/don't have a relationship with the person. e.g.: 'because of ignorance', 'we are not so close', 
                'signal' - Observing (or not) a bracelet, ink on my thumb, or a calendar. e.g.: 'he saw my ink', 'he didn't see my bracelet'
                'type' - They know I'm a good person that cares about my health and community (or a bad person that doesn't). e.g.: 'I always attend such activities', 'knows i cannot miss', 'because she knows i love such things'
                'circumstances' - There were circumstances that prevented me from getting dewormed. e.g.: 'he knows i didnt get information', 'they know am blind' 
                'other' - None of the above fit well.
                """ 
                }, 
            {"role": "user", "content": "The string to classify is: " + reason}, 

        ]
    )
    print([reason, completion.choices[0].message.content])

    # write to csv and append:
    if write_csv:
        with open("temp-data/sob-other-reasons-gpt-3.5-turbo.csv", "a") as f:
            f.write(reason + "," + completion.choices[0].message.content)

    return completion.choices[0].message.content

# classify_reason_with_gpt("I told them I was dewormed")
# classify_reason_with_gpt("a family member")
# classify_reason_with_gpt("according to how she knows me")
# classify_reason_with_gpt("all went")
# classify_reason_with_gpt("almost eveyone comes")

sob_reason_csv = pd.read_csv("temp-data/sob-other-reasons.csv")

sob_reason_csv["gpt-3.5-turbo"] = sob_reason_csv["second.order.reason"].apply(classify_reason_with_gpt, write_csv = True)

sob_reason_csv.to_csv("temp-data/sob-other-reasons-gpt-3.5-turbo-single-pass.csv", index=False)

